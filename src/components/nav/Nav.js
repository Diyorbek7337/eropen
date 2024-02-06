import React from 'react';
import Logo from '../../image/logo.png';
import { FaGlobe} from 'react-icons/fa'
import  {useState, useEffect}  from 'react';
import { useTranslation } from 'react-i18next';

const Nav = ({changeLan}) => {
    const [open, setOpen] = useState(true)
    const [lang, setLang] = useState(true)
    const [ show, setShow] = useState(false)
    useEffect(() => {
      setShow(true)
  },[])
    const { t } = useTranslation()
    const lockMobile = () => {
    if (open)
      document.body.style.overflow = 'unset';
     else 
     document.body.style.overflow = 'hidden';
    }
    lockMobile()

    return (
        <>
            <header id="header" className={show ? "header --show-item --show" : "header lock-padding --show-item"}>
          <div className="header__row">
            <img src={Logo} alt="Dora Line" className='logo' />
            <div className="header__menu menu">
              <nav className={open ? "menu__navbar" : "menu__navbar --active"}>
                <ul className="menu__list">
                  <li className={open ? "menu__list-item menu__list-item-1" : "menu__list-item menu__list-item-1 --active"}>
                    <a className="menu__link" href="#mains">{t("nav.home")}</a>
                  </li>
                  <li className={open ? "menu__list-item menu__list-item-2" : "menu__list-item menu__list-item-2 --active"}>
                    <a className="menu__link" href="#contact">{t("nav.connect")}</a>
                  </li>
                  <li className={open ? "menu__list-item menu__list-item-3" : "menu__list-item menu__list-item-3 --active"}>
                    <a className="menu__link popup-link" href="#bat">{t("nav.more")}</a>
                  </li>
                  <li className={open ? "menu__list-item menu__list-item-4" : "menu__list-item menu__list-item-4 --active"}>
                    <div className="menu__lang-list--row">
                      <a className="menu__link menu__link-4" href='#'  onClick={() => setLang(!lang)}><FaGlobe/> {t("nav.language")}</a>
                      {/* <a className="menu__link menu__link-4" href='#' onClick={() => setLang(!lang)}><{t("nav.language")}</a> */}
                      <ul className={lang ? "menu__lang-list lang-list" : "menu__lang-list lang-list --active"}>
                        <li className="lang-list__item">
                          <a href="#" className="lang-list__link --active"  onClick={() => changeLan('Uz')}>O'zbekcha</a>
                        </li>
                        <li className="lang-list__item">
                          <a href="#" className="lang-list__link" onClick={() => changeLan('Ru')}>Русский</a>
                        </li>
                        <li className="lang-list__item">
                          <a href="#" className="lang-list__link" onClick={() => changeLan('En')}>English</a>
                        </li>
                      </ul>
                    </div>
                  </li>
                </ul>
              </nav>
              <div className={open ? "menu__burger-menu burger-menu" : "menu__burger-menu burger-menu --active"} onClick={() => setOpen(!open)}>
                <span className="burger-menu__line"></span>
              </div>
          </div>
          </div>
        </header>
        </>
    );
};

export default Nav;