import React from 'react';
import Logo from '../../image/logo.png'
import './loading.css'

const Loading = ({open}) => {
    const isOpen = open ? 'loader active' : 'loader'
    return (
        
            <div className='wrapper'>
                <div className={isOpen}>
                    <img src={Logo} alt="Dora Line" className='loader-image'/>
                </div>
            </div>
        
    );
};

export default Loading;